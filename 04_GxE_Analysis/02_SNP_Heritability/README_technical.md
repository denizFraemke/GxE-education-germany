# 02_SNP_Heritability — Analysis S5 (GREML): technical / reproducibility reference

> **This is the full technical reference** — pipeline, workflow commands,
> methods detail, SLURM steps, diagnostics and output tree. For the short
> summary (analysis hierarchy and limitations) see
> [`README.md`](README.md). Departures from the registered plan are in
> [`../../Plan_deviations.md`](../../Plan_deviations.md) (§7a–§7i).
> The estimates themselves are reported in the manuscript, not here.

SNP-based heritability of years of education by **Region × Reunification**
(East/West × pre/post-1990) in the harmonized pooled sample, per Analysis
S5 of `ANALYSIS_PLAN.md`. Three additional secondary / exploratory
comparison tracks (R1 Region, G1 recorded gender/sex, RG1 Region ×
gender/sex) run alongside the primary on the same source analytic
sample before track-specific kinship pruning.

The pipeline takes the Phase-1 harmonized cohort bfiles from
`01_Genotype/Harmonized/` (re-filtered to imputation R² ≥ 0.9), merges
them into a single pooled bfile, and produces four GREML
specifications run on the same source analytic sample (retained N
differs across specifications after track-specific kinship pruning):

- **Primary cell-level analysis:** *block-diagonal GREML with cohort
  fixed effects* (steps `03d/04d/05d/06d`). Per-cohort GRMs are
  combined block-diagonally so cross-cohort GRM entries are zero by
  design; V_G is estimated from within-cohort genomic relatedness only
  while cohort-level mean differences in years of education are
  adjusted in `--covar` (`gender_int + cohort_int`). See
  `Plan_deviations.md` §7g for the diagnostic that established this
  specification.

- **Sensitivity 1** (`03/04/05`): single pooled cross-cohort GRM per
  cell at a liberal kinship filter `--grm-cutoff 0.20`. This filter
  can retain substantial relatedness and is not equivalent to a
  conventional unrelated-sample GREML; it is reported as a sensitivity
  rather than the primary.

- **Sensitivity 2** (`03c/04c/05c/06c`): per-cohort GREML using
  cohort-specific GRMs at `--grm-cutoff 0.05`. Answers "is any single
  cohort driving the pooled estimate?".

- **Secondary / exploratory comparison tracks** (`03e/04d/05d/06d` via
  the three `run_blockdiag_*.sh` wrappers, dispatched as `--step xtra`):
  same block-diagonal design as the cell-level primary, but applied to
  three different stratifications — R1 Region (East/West, collapsed
  across pre/post-1990), G1 recorded gender/sex, and RG1 Region ×
  recorded gender/sex. Uses a *subset-then-prune* GRM extractor (03e)
  that avoids removing individuals solely because they are related to
  someone who falls into a different comparison stratum. The validated
  cell-level primary is unchanged and protected by
  `regression_check_blockdiag_primary.sh`, which aborts the secondary
  tracks if the primary's h²/SE drifts by > 0.001. See
  `Plan_deviations.md` §7h.

Cross-cell h² heterogeneity is tested by Cochran's Q + pairwise Z on
the four per-cell h² estimates from the primary specification.
Precision is summarized by the GCTA REML SEs from the `.hsq` output;
analytical Visscher (2014) power estimates are computed Tardis-side as
a diagnostic but are not reported in the headline (the formula's
unrelated-sample assumption is violated by the sensitivity-1 cutoff and
by retained within-cohort relatedness — see *Methods → Power calc*).

> **How to read this section.** The primary cell-level analysis is the
> block-diagonal GREML with cohort fixed effects (steps
> `03d/04d/05d/06d`; `Plan_deviations.md` §7g). It estimates SNP-based
> variance from within-cohort genomic relatedness only, while adjusting
> for cohort-level mean differences in years of education. The
> pooled-GRM cell-level analysis (`GRM_CUTOFF=0.20`) and the per-cohort
> GREML are reported as sensitivities.
>
> Three secondary / exploratory comparison tracks (R1 Region, G1
> recorded gender/sex, RG1 Region × gender/sex) run alongside the
> cell-level primary as supplementary checks. See `Plan_deviations.md`
> §7h and `output/final/h2_{region_main,gender_main,region_gender}_summary.md`.

> **Execution.** This section is **uploaded to Tardis and run there**.
> The Mac side only writes the phenotype/stratum TSVs and pushes them +
> the code over rsync. There is no Mac-side compute beyond
> `export_stratum_phenotype.R`. After the pipeline finishes, the
> Mac-side `summarize_results.R` builds the final report + plots from
> the Tardis-produced TSVs (no genotype data needed). See *Workflow*
> below.

---

## Build provenance and outputs

**Build provenance.** GCTA 1.95.1; `GRM_CUTOFF=0.20` (pooled-GRM sensitivity), `GRM_CUTOFF_PER_COHORT=0.05` (per-cohort sensitivity and block-diagonal primary; applied within each cohort before block stacking); `R²_THRESHOLD=0.9` with the SOEP chr1 fallback (see *TL;DR → Variant set*); analytic sample from `Combined_Harmonized_20260513.rds`; covariates `crossPC1..20 + birth_year_c + birth_year_c_sq` in `--qcovar` and `gender_int` in `--covar` (plus `cohort_int` in `--covar` for the block-diagonal primary track; see §7g).

Block-diagonal primary outputs: `output/final/h2_blockdiag.tsv`, `h2_blockdiag_heterogeneity.tsv`, `h2_blockdiag_summary.md`. Pooled-GRM sensitivity outputs: `output/final/h2_snp.tsv`, `heterogeneity.tsv`. Per-cohort sensitivity outputs: `output/final/h2_per_cohort.tsv`, `h2_per_cohort_summary.md`. Secondary-track outputs (one set per track, TRACK ∈ {`region_main`, `gender_main`, `region_gender`}): `output/final/h2_<TRACK>.tsv`, `h2_<TRACK>_heterogeneity.tsv`, `h2_<TRACK>_summary.md`. Forest plots (Mac-side, produced by `summarize_results.R`): **one single combined multi-page PDF** `output/final/h2_forest_plots.pdf` (one analysis per page, in order: primary block-diagonal → pooled-GRM sensitivity → per-cohort sensitivity → R1 → G1 → RG1) plus per-track PNGs (`h2_blockdiag_forest_plot.png`, `h2_forest_plot.png`, `h2_per_cohort_forest_plot.png`, `h2_<TRACK>_forest_plot.png`) useful for embedding individual plots into reports / slides. Diagnostic outputs: `output/final/h2_blockdiag_diagnostic.{tsv,md}` (single-block replication + no-cohort-covar comparison). All per-stratum TSVs (cell-level primary, both sensitivities, and the three secondary tracks) carry the underlying variance components — `V_G`, `V_G_SE`, `V_e`, `V_e_SE`, `Vp` — in columns 5–9; the per-track `summary.md` files include a Variance-components table parallel to the h² table (per `ANALYSIS_PLAN.md` §S5 and `Plan_deviations.md` §7i).

---

## TL;DR

- **4-cell stratification:** East/West × pre/post-1990 reunification (cutoff = turned 15 in/after 1990 ⇔ BY ≥ 1975).
- **Variant set:** cross-cohort intersection at imputation R² ≥ 0.9, MAF ≥ 0.01, autosomes only (HWE filter skipped — see `Plan_deviations.md` §7d). **Exception:** SOEP chr1 INFO file was never delivered by the data provider; the upstream `03a_filter_imputation_quality.sh` falls back to including all chr1 variants from the SOEP harmonized `.bim` at the mildQC R² > 0.1 threshold. Reflected in `manifest.tsv` as `soep_chr1_info=fallback_mildQC_R2_0.1`.
- **SNP set, by track.**
  - **Primary block-diagonal track:** each (cohort × cell) block inherits the cohort-level SNP set and within-cohort allele-frequency standardization from the per-cohort GRMs (step 03c), before cohort × cell subsetting (step 03d).
  - **Pooled-GRM sensitivity (cutoff 0.20):** GCTA's `--maf 0.01` filter is applied within each cell's `--keep` subset, so each cell's pooled GRM uses SNPs at MAF ≥ 0.01 in that cell. The four cell GRMs may therefore use slightly different SNP sets. A sensitivity using a fixed pooled-MAF global SNP list is not implemented.
- **Relatedness handling, by track.**
  - **Primary block-diagonal track:** within each cohort, `--grm-cutoff 0.05` is applied (the published convention). Cross-cohort GRM entries are zero by construction.
  - **Pooled-GRM sensitivity:** `--grm-cutoff 0.20`, a liberal kinship filter that can retain substantial relatedness and is not equivalent to a conventional unrelated-sample GREML. The plan calls for 0.025; on the pooled cross-cohort GRM, thresholds in the 0.025–0.10 range remove the minority cohort of each Region × Reunification cell almost entirely, so 0.20 is used instead (see `Plan_deviations.md` §7e). The kinship-pruned GRMs carry the filename suffix `_unrel` regardless of the cutoff actually applied.
- **Covariates, by track.**
  - **Numeric (`--qcovar`, all tracks):** `crossPC1..20`, `birth_year_c`, `birth_year_c_sq` (the BY adjustment is in `Plan_deviations.md` §7f).
  - **Categorical (`--covar`).** Primary block-diagonal track: `gender_int + cohort_int` from `gender_cohort.covar`. `cohort_int` is categorical (1=BASE-II, 2=SHIP, 3=SOEP, 4=TwinLife) and adjusts cohort-level mean differences in years of education that would otherwise load into V_e under the block-diagonal design (cross-cohort A entries are zero, so cohort means cannot be absorbed by V_G; see §7g). Pooled-GRM sensitivity and per-cohort sensitivity use the `gender_int`-only specification from `gender.covar`.
  - **Secondary tracks** (`Plan_deviations.md` §7h). R1 Region: same `gender_cohort.covar` as the primary. G1 recorded gender/sex: `region_cohort.covar` (region_int + cohort_int — gender is the stratifier, so it's not also in `--covar`). RG1 Region × gender/sex: `cohort_only.covar` (cohort_int only — both stratifying factors define the stratum).
- **Phenotype:** `education` in raw years (V_G interpretable in years²).
- **Cross-cell h² test:** **Cochran's Q heterogeneity test** on the K cell h² estimates from the primary specification + pairwise Z-tests. The `gcta --reml --mgrm` joint LRT specified in the plan segfaults on disjoint-sample inputs — see `Plan_deviations.md` §7b. Step 07 implements the Q-based test.
- **Power:** **not reported analytically.** The Visscher (2014) `SE(h²)` formula assumes an unrelated sample; at the pooled-GRM-sensitivity cutoff of 0.20 it materially underestimates the SE in cells that retain meaningful relatedness. Precision is summarized by the GCTA REML SE from the `.hsq` output. The Visscher script runs Tardis-side as a diagnostic but its output (`power.tsv`) is not surfaced in the report.
- **Per-cohort sensitivity track:** steps `03c/04c/05c/06c` fit GREML per cohort using cohort-specific GRMs at `--grm-cutoff 0.05`. Answers "is one cohort driving the result?".
- **Secondary / exploratory comparison tracks (R1, G1, RG1):** same block-diagonal design as the cell-level primary, but applied to three other stratifications (R1 Region, G1 recorded gender/sex, RG1 Region × gender/sex) using a *subset-then-prune* GRM extractor (`03e_extract_cohort_x_stratum_grms.{sh,slurm}`) — the per-(cohort × stratum) sub-GRM is subset from the pre-cutoff per-cohort GRM, then `--grm-cutoff 0.05` is applied *within* the block. This avoids removing individuals solely because they are related to someone who falls into a different comparison stratum. Guarded by `scripts/regression_check_blockdiag_primary.sh` (aborts the secondary tracks if the primary's h²/SE drifts by > 0.001). Dispatched as `bash run_pipeline.sh --step xtra`; see `Plan_deviations.md` §7h. **Treat results as exploratory / sensitivity, not as standalone tests of group differences in h².**

See `Plan_deviations.md` (repo root) for departures from the literal plan.

---

## Workflow — Mac side, then Tardis

### One-time setup (Mac)

```sh
# Tardis SSH key — verify there's no password prompt:
ssh ${CLUSTER_USER}@${CLUSTER_HOST} echo ok
# If it prompts: ssh-copy-id ${CLUSTER_USER}@${CLUSTER_HOST}
```

`user_config.sh` is generated interactively on first run of any
`transfer_*.sh` script (or by `cp user_config.sh.example user_config.sh`).
It is **gitignored**.

### Run

| Step | Where | Command | What it does |
|---|---|---|---|
| 1 | Mac | `Rscript 04_GxE_Analysis/02_SNP_Heritability/export_stratum_phenotype.R` | Reads latest `Combined_Harmonized_<YYYYMMDD>.rds`; writes `education.phen`, `crosspcs.qcovar`, `gender.covar`, `gender_cohort.covar`, `region_cohort.covar`, `cohort_only.covar`, `keep/{<cell-stratum>,region_*,gender_*,region_gender_*}.iids`, and `manifest.tsv` to `${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/data/`. Aggregate-only stdout. |
| 2 | Mac | `bash 04_GxE_Analysis/02_SNP_Heritability/transfer_phenotype_to_tardis.sh` | rsync the `data/` directory above → `${TARDIS_PROJECT_ROOT}/data/`. |
| 3 | Mac | `bash 04_GxE_Analysis/02_SNP_Heritability/transfer_pipeline_to_tardis.sh` | rsync the **code** → `${TARDIS_PROJECT_ROOT}/`. |
| 4 | Tardis | `cd ${TARDIS_PROJECT_ROOT} && bash run_pipeline.sh` | Runs steps 00 → 09, then the per-cohort sensitivity (03c → 06c), the block-diagonal primary cell-level track (03d → 06d), and the secondary / exploratory comparison tracks (`xtra`: regression check + R1 + G1 + RG1). Each step is idempotent. |
| 5 | Mac | `bash 04_GxE_Analysis/02_SNP_Heritability/transfer_results_from_tardis.sh` | rsync the final TSVs back to the shared volume. |
| 6 | Mac | `Rscript 04_GxE_Analysis/02_SNP_Heritability/summarize_results.R` | Builds `summary_report.md`, **one combined multi-page PDF `h2_forest_plots.pdf`** with one forest plot per analysis track (primary → pooled-GRM sensitivity → per-cohort → R1 → G1 → RG1), plus per-track PNGs for embedding (`h2_forest_plot.png`, `h2_blockdiag_forest_plot.png`, `h2_per_cohort_forest_plot.png`, `h2_region_main_forest_plot.png`, `h2_gender_main_forest_plot.png`, `h2_region_gender_forest_plot.png`). **Does not need any SNP data.** |

After step 6 completes you have a self-contained results bundle on the
shared volume — safe to delete cohort-level SNP data on Tardis.

### Resume / single step (on Tardis)

```sh
bash run_pipeline.sh --step 05
bash run_pipeline.sh --from 03 --to 06
bash run_pipeline.sh --step xtra     # regression check + R1 + G1 + RG1
```

---

## Pipeline overview

Numbered steps under `scripts/` all run on Tardis (sourced from `config.sh`).

| Step | Script | Purpose |
|------|--------|---------|
| 00 | `00_check_inputs.sh` | Verify gcta64/plink2/plink1.9, upstream Phase-1 bfiles + info files, phenotype + covariates + per-stratum keep files. |
| 01 | `01_rebuild_bfiles_info090.sh` | Re-run upstream `01_Genotype/Harmonized/scripts/03a_filter_imputation_quality.sh` at `R2_THRESHOLD_STRICT=0.9`; cross-cohort intersection; `plink2 --extract` per cohort from the Phase-1 unrestricted bfiles → `{COHORT}_harmonized_r2_0.9.{bed,bim,fam}`. |
| 02 | `02_pool_cohorts.sh` | `plink 1.9 --merge-list` the 5 INFO ≥ .90 bfiles → `pooled.{bed,bim,fam}` with missnp-retry. |
| 03 | `03_build_stratum_grms.sh` + `03_build_stratum_grms.slurm` | SLURM array, one task per stratum: `gcta64 --bfile pooled --keep <stratum>.gcta.iids --autosome --maf 0.01 --make-grm-bin`. Compute-node defaults: `partition=short`, `--exclusive`, `--cpus-per-task=16`, `--mem=0`, `--time=8:00:00`. |
| 04 | `04_unrelated_filter.sh` + `04_unrelated_filter.slurm` | Per stratum: `gcta64 --grm <stratum> --grm-cutoff 0.20 --make-grm-bin --out grm_<stratum>_unrel`. Despite the `unrelated_filter` name, the operation at cutoff 0.20 is a liberal kinship filter that can retain substantial relatedness, not a conventional unrelated-sample filter. Output is consumed by the pooled-GRM sensitivity (step 05). |
| 05 | `05_fit_greml_main.sh` + `05_fit_greml_main.slurm` | Submits 8 per-stratum standalone REML fits (4 cells × constrained / `--reml-no-constrain`) as a SLURM array. |
| 07 | `07_lrt.R` | Reads per-cell standalone `.hsq` files; computes **Cochran's Q + pairwise Z + I²** → `output/reml/heterogeneity.tsv` (compatibility symlink `lrt.tsv`). |
| 08 | `08_power_calc.R` | Per stratum: empirical Var(off-diag-GRM) from each `*_unrel.grm.bin`; Visscher 2014 SE(h²) and chi² LRT power across the h² grid → `output/reml/power.tsv`. |
| 09 | `09_assemble_output.R` | Builds `output/final/{h2_snp,heterogeneity,lrt,power,manifest,dashboard_row}.{tsv,md}` from the per-step outputs. |

(There is no step 06; the 05 → 07 numbering gap is intentional.)

### Per-cohort sensitivity track

Parallel to the cell-level pipeline (steps 03/04/05), four additional steps fit GREML *per cohort* instead of per stratum. Motivation: (a) the cell-level pooled GRM carries cross-cohort fingerprint that cutoff-tuning cannot eliminate (see §7e); (b) "is any single cohort driving the headline result?" is best answered with per-cohort h². Per-cohort GRMs eliminate both the fingerprint problem and the minority-cohort-wipe-out problem by construction.

| Step | Script | Purpose |
|------|--------|---------|
| 03c | `03c_build_per_cohort_grms.{sh,slurm}` | Build one GRM per cohort (BASE-II, SHIP = SHIP-0 + SHIP-Td, SOEP, TwinLife) from pooled.bfile with `--keep <cohort-IIDs>`. Within-keep MAF standardization → no cross-cohort fingerprint. |
| 04c | `04c_unrelated_filter_per_cohort.{sh,slurm}` | `gcta --grm-cutoff 0.05` per cohort (published convention; no cross-cohort fingerprint inside one cohort). |
| 05c | `05c_fit_greml_per_cohort.{sh,slurm}` | One REML per cohort (full cohort), plus cohort × region and cohort × stratum sub-fits via `--keep` for any sub-sample with N ≥ 100. |
| 06c | `06c_assemble_per_cohort.R` | `output/final/h2_per_cohort.tsv` + `h2_per_cohort_summary.md`. |

### Block-diagonal cell-level track (primary cell-level estimator; with cohort fixed effects)

| Step | Script | Purpose |
|------|--------|---------|
| 03d | `03d_extract_cohort_x_cell_grms.{sh,slurm}` | For each (cohort × cell) sub-sample with N ≥ `MIN_BLOCKCOMP_N` (= 30), subset the per-cohort unrel GRM via `gcta --grm-bin --keep --make-grm-bin`. Output: 11–13 `grm_blockcomp_<COHORT>_<CELL>.grm.{bin,id,N.bin}`. Cheap (no genotype I/O). |
| 04d | `04d_build_block_diagonal_cell_grms.sh` | Per cell: combine the cohort blocks into one block-diagonal GRM via `Rscript R/build_block_diagonal_grm.R`. Output: 4 `grm_blockdiag_<CELL>.grm.{bin,id,N.bin}`. Runs on the login node (no SLURM). |
| 05d | `05d_fit_greml_block_diagonal.{sh,slurm}` | 8 REML fits (4 cells × constrained / `--reml-no-constrain`) on the block-diagonal GRMs, using `gender_cohort.covar` in `--covar` (gender + cohort indicator). Reuses the FID-alignment block from step 05; aligns the block-diagonal covar file as well. |
| 06d | `06d_assemble_block_diagonal.R` | `output/final/h2_blockdiag.tsv` + `h2_blockdiag_heterogeneity.tsv` + `h2_blockdiag_summary.md`. |

These steps are idempotent and ordered alphanumerically after step 06c in `run_pipeline.sh`. Run the per-cohort + block-diagonal sensitivity tracks with `bash run_pipeline.sh --from 03c --to 06d` after the main pipeline has finished (or just `bash run_pipeline.sh` to run everything).

### Secondary comparison tracks (R1 / G1 / RG1; exploratory; `--step xtra`; §7h)

Three additional block-diagonal GREML comparisons run on the same source analytic sample (before track-specific kinship pruning), dispatched together by a single new pipeline step (`xtra` in `run_pipeline.sh`):

| Track | Strata | `--covar` | Wrapper |
|---|---|---|---|
| R1 | `region_east`, `region_west` | `gender_cohort.covar` | `scripts/run_blockdiag_region_main.sh` |
| G1 | `gender_male`, `gender_female` | `region_cohort.covar` | `scripts/run_blockdiag_gender_main.sh` |
| RG1 | `region_gender_{east,west}_{male,female}` | `cohort_only.covar` | `scripts/run_blockdiag_region_gender.sh` |

Each wrapper sets `TRACK`, `STRATA_LIST_NAMES`, `BLOCKDIAG_OUT_SUBDIR`, and `BLOCKDIAG_COVAR_FILE_OVERRIDE` via env vars, then dispatches the chain `03e → 04d → 05d → 06d` with those settings. The cell-level primary chain (`03d → 04d → 05d → 06d`) is unaffected — `04d/05d/06d` default to the validated primary stratum list and covar file when those env vars are unset.

Before any secondary track runs, `scripts/regression_check_blockdiag_primary.sh` reads `output/final/h2_blockdiag.tsv` and verifies the four primary cells' constrained-fit h² and REML SE match the validated reference values (East × pre-1990 0.282 / 0.050; East × post-1990 0.040 / 0.428; West × pre-1990 0.281 / 0.172; West × post-1990 0.000 / 0.347) within ±0.001. If the check fails, the secondary tracks do **not** run, and the suspected upstream change is investigated first.

| Step | Script | Purpose |
|---|---|---|
| xtra-regcheck | `scripts/regression_check_blockdiag_primary.sh` | Guard: primary cell-level h²/SE unchanged within ±0.001. |
| 03e | `03e_extract_cohort_x_stratum_grms.{sh,slurm}` | For each (cohort × stratum) with pre-cutoff N ≥ `MIN_BLOCKCOMP_N` (= 30): subset the pre-cutoff per-cohort GRM (`grm_cohort_<COHORT>` from step 03c) using an FID-aligned `(cohort × stratum)` keep file, then `--grm-cutoff 0.05` *within* that sub-block. Output: `grm_blockcomp_<COHORT>_<STRATUM>.{bin,id,N.bin}`. The FID-aligned keep is built per task by intersecting the cohort's `.grm.id` with the stratum's `.iids` on IID alone — GCTA `--keep` matches on (FID, IID), so the stratum keep's source-data FIDs must be rewritten to the cohort-pooled FIDs first. |
| 04d | `04d_build_block_diagonal_cell_grms.sh` | Same script as the primary, but `STRATA_LIST_NAMES` is set to the secondary-track stratum list. Combines per-(cohort × stratum) blocks into one `grm_blockdiag_<STRATUM>` per stratum. |
| 05d | `05d_fit_greml_block_diagonal.{sh,slurm}` | Same script as the primary, but with `BLOCKDIAG_OUT_SUBDIR` (subdir under `output/reml/`) and `BLOCKDIAG_COVAR_FILE_OVERRIDE` (the track's covar file) set by the wrapper. |
| 06d | `06d_assemble_block_diagonal.R` | Same script as the primary, parameterised by `TRACK` (output filename suffix), `STRATA_LIST_NAMES`, `BLOCKDIAG_OUT_SUBDIR`, `BLOCKDIAG_COVAR_FILE_OVERRIDE`. Emits per-track `h2_<TRACK>.{tsv,md}`, `h2_<TRACK>_heterogeneity.tsv` plus the per-track diagnostic blocks (cohort composition, single-cohort-dominance flag at >80%, boundary flag, edu mean/SD/Var, mean centered birth year). |

**Why subset-then-prune?** The cell-level primary uses prune-then-subset (per-cohort GRM filtered cohort-wide at `--grm-cutoff 0.05` in step 04c, then subset to cells in step 03d). For tracks like G1, that order can remove an individual at the cohort-wide step solely because of a within-cohort relative who would have ended up in a different stratum (e.g. a cross-sex sib pair in TwinLife under G1). Subset-then-prune avoids removing individuals solely because they are related to someone who falls into a different comparison stratum; the within-block `--grm-cutoff 0.05` still applies, so each retained block remains unrelated within itself. The trade-off is a more aggressive within-block prune, quantified empirically: BASE-II loses 0%, TwinLife 23–30%. The cell-level primary deliberately stays on the prune-then-subset order, which is what its validated h² values are based on; see `Plan_deviations.md` §7h for the full design discussion.

After the Tardis side finishes, **on the Mac:**

| Helper | What it does |
|---|---|
| `transfer_results_from_tardis.sh` | rsync `${TARDIS_PROJECT_ROOT}/output/final/` to `${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/output/final/`. |
| `summarize_results.R` | Builds `summary_report.md`, **one combined multi-page PDF `h2_forest_plots.pdf`** containing one forest plot per analysis track (six pages total), plus per-track PNGs for embedding, and validates the bundle. **Reads only TSVs — no SNP data needed.** Pre-deletion sanity check. |

---

## Methods

### Stratification

`Region` from the merged-RDS column `east_west` (character, `"east"/"west"`).
`Time` from `birth_year` via the "turned-15-in-or-after-1990" rule
(ANALYSIS_PLAN.md §German reunification and A2/B2): `BY < 1975 → pre1990`,
`BY ≥ 1975 → post1990`. There is no pre-computed reunification indicator
in the merged RDS — `export_stratum_phenotype.R` derives it from
`birth_year`.

The secondary comparison tracks (§7h) define additional strata on the
same source analytic sample (before track-specific kinship pruning):

- **R1 — Region main:** `region_east`, `region_west`, derived from the
  same `east_west` column but collapsed across pre/post-1990. Note this
  is *not* the union of the two cell-level Region × Reunification
  strata as far as the GRM goes; the GRMs are rebuilt by 03e
  independently.
- **G1 — recorded gender/sex main:** `gender_male`, `gender_female`,
  derived from a harmonised `gender_norm` column
  (`Male / Female`-only after normalisation; rows with non-standard
  values are dropped from G1 and RG1).
- **RG1 — Region × recorded gender/sex:** four strata
  (`region_gender_east_male`, `region_gender_west_male`,
  `region_gender_east_female`, `region_gender_west_female`), defined as
  the intersection of `region_*` and `gender_*` above.

Stable categorical codings used in the new covar files:
`region_int` (1 = East, 2 = West), `gender_int` (1 = Male, 2 = Female),
`cohort_int` (1 = BASE-II, 2 = SHIP, 3 = SOEP, 4 = TwinLife).
"Recorded gender/sex" is the harmonised label from the source cohorts;
it is not a measure of gender identity or biological sex per se, and
the comparison should be read accordingly (see §7h).

### Variant set

INFO ≥ 0.9 is enforced via the upstream pipeline's pass-list mechanism
(`R2_THRESHOLD_STRICT=0.9`), intersected across the 5 cohorts. MAF ≥
0.01 is applied by `gcta --make-grm-bin`. HWE filter is **not**
applied — see `Plan_deviations.md` §7d. The **SOEP chr1 INFO file was
never delivered by the data provider**; the upstream
`03a_filter_imputation_quality.sh` falls back to including all chr1
variants from the SOEP harmonized `.bim` at the mildQC R² > 0.1
threshold. Reflected in `manifest.tsv` as
`soep_chr1_info=fallback_mildQC_R2_0.1`. The MAF filter applies per
track: the pooled-GRM sensitivity computes MAF within each cell's
`--keep` subset; the primary block-diagonal track inherits the
within-cohort MAF / standardization from the per-cohort GRMs (step
03c) before cohort × cell subsetting.

### Covariates

Numeric (`--qcovar`, all tracks): `crossPC1..20`, `birth_year_c`,
`birth_year_c_sq`, where
`birth_year_c = birth_year − mean(birth_year)` on the analytic sample
(centering is for numerical conditioning of the quadratic term;
see `Plan_deviations.md` §7f).

Categorical (`--covar`), **primary block-diagonal track**:
`gender_int + cohort_int` from `gender_cohort.covar`. `cohort_int` is
categorical (1=BASE-II, 2=SHIP, 3=SOEP, 4=TwinLife) and is included
specifically to adjust cohort-level mean differences in years of
education that would otherwise load into V_e under the block-diagonal
design — see `Plan_deviations.md` §7g.

Categorical (`--covar`), **pooled-GRM sensitivity and per-cohort
sensitivity**: `gender_int`-only specification from `gender.covar`
(1=Male, 2=Female).

Categorical (`--covar`), **secondary comparison tracks (§7h)**: the
covar file is chosen per track so that the stratifying factor itself
is *not* also in `--covar` (it would be collinear with the
block-diagonal stratum structure):

- **R1 Region main** uses `gender_cohort.covar` — same file as the
  cell-level primary (`gender_int + cohort_int`).
- **G1 recorded gender/sex main** uses `region_cohort.covar`
  (`region_int + cohort_int`). Gender is the stratifier, so it isn't in
  `--covar`.
- **RG1 Region × recorded gender/sex** uses `cohort_only.covar`
  (`cohort_int` only). Both Region and Gender define the stratum, so
  neither is in `--covar`.

All three retain the cohort indicator. The block-diagonal GRMs have
cross-cohort entries set to zero by construction, so cohort-level mean
differences in years of education cannot be absorbed into V_G and must
be modelled as fixed effects to keep V_e from inflating — same logic
as §7g.

### Per-cell REML

**Primary: block-diagonal GREML with cohort fixed effects.** For each
cell k, build the cell GRM as a block-diagonal matrix of per-(cohort ×
cell) sub-GRMs:

```
A_k = blockdiag( A_{k,BASEII}, A_{k,SHIP}, A_{k,SOEP}, A_{k,TwinLife} )
```

with within-cohort allele-frequency standardization and zero
cross-cohort entries by construction. Fit

```
y_k = X_k β_k + g_k + e_k
g_k ~ N(0, σ²_G_k · A_k)
e_k ~ N(0, σ²_e_k · I_k)
```

via `gcta --reml --grm-bin grm_blockdiag_<k> --covar
gender_cohort.covar --qcovar crosspcs.qcovar`. X_k therefore includes
cohort fixed effects (dummy-coded internally by GCTA), partialling
out cohort-level mean differences in years of education before
variance partitioning. V_G is identified from within-cohort genomic
relatedness signal only. h²_SNP_k = σ²_G_k / (σ²_G_k + σ²_e_k); SE
from GCTA's AI-REML output.

**Sensitivity 1: single pooled cross-cohort GRM per cell.** Same model
form as above, but the cell GRM `A_k` is the single pooled GRM built
from the merged bfile with `gcta --bfile pooled --keep <stratum>.iids
--make-grm-bin` and then filtered with `--grm-cutoff 0.20`. Cross-cohort
entries are non-zero in this design, so cohort-mean differences are
partly absorbed into V_G via inflated relatedness; a cohort indicator
is therefore **not** added to `--covar` here, and the gender-only
specification keeps this track comparable across cutoffs.
The 0.20 cutoff is a liberal kinship filter
that can retain substantial relatedness; this specification is
reported as a sensitivity, not the primary.

**Sensitivity 2: per-cohort GREML.** One REML per cohort using
cohort-specific GRMs at `--grm-cutoff 0.05`; `--covar` retains the
gender-only specification (cohort is constant within a fit).

**Secondary tracks (R1 / G1 / RG1): subset-then-prune block extraction.**
Same model form as the primary, but the per-(cohort × stratum) sub-GRMs
that make up each `A_k` are extracted by a different script (03e) and a
different ordering. For each (cohort × stratum) combination with
pre-cutoff intersection N ≥ `MIN_BLOCKCOMP_N` (= 30), 03e runs two GCTA
calls in sequence:

1. `gcta --grm-bin grm_cohort_<COHORT> --keep <FID-aligned (cohort × stratum) keep> --make-grm-bin --out <intermediate>`
   — subset the pre-cutoff per-cohort GRM (`grm_cohort_<COHORT>` from
   step 03c, *not* the post-cutoff `_unrel` version) to the
   (cohort × stratum) sub-sample.
2. `gcta --grm-bin <intermediate> --grm-cutoff 0.05 --make-grm-bin --out grm_blockcomp_<COHORT>_<STRATUM>`
   — apply the conventional 0.05 within-block unrelated filter.

The intermediate GRM is cleaned up by a trap in the SLURM worker. The
FID-aligned keep file in step (1) is built by intersecting the cohort's
`.grm.id` (column 1 = canonical FID, column 2 = IID) with the stratum
keep file on IID alone, then emitting `(FID_from_cohort_GRM, IID)`
pairs — GCTA's `--keep` matches on `(FID, IID)` and the stratum keep
file inherits source-data FIDs that don't match the cohort-pooled
`pooled.fam` FIDs for BASE-II and TwinLife.

This *subset-then-prune* order differs from the cell-level primary's
*prune-then-subset* (the primary uses the post-cutoff cohort GRM
`grm_cohort_<COHORT>_unrel` from step 04c, then subsets to cells in
step 03d, with no additional within-block cutoff). The secondary
tracks deliberately adopt the inverse order so that an individual is
not removed at the cohort-wide step solely because of a within-cohort
relative who would have ended up in a different comparison stratum
(e.g. a cross-sex sib pair in TwinLife under G1). The within-block
`--grm-cutoff 0.05` still applies, so each retained block remains
unrelated within itself. The primary stays on prune-then-subset, which
is what its validated h² values are based on. See
`Plan_deviations.md` §7h for the full rationale and `scripts/regression_check_blockdiag_primary.sh` for
the guard that enforces it.

The combined block-diagonal GRM for each secondary-track stratum is
then built by `04d` (same script as the primary, parameterised by
`STRATA_LIST_NAMES`) and fit by `05d` (same script, with
`BLOCKDIAG_OUT_SUBDIR` and `BLOCKDIAG_COVAR_FILE_OVERRIDE` set per
track). `06d` assembles per-track outputs (`h2_<TRACK>.{tsv,md}`,
`h2_<TRACK>_heterogeneity.tsv`) plus per-track diagnostics: cohort
composition table, single-cohort-dominance flag (>80%), boundary flag
on near-0 / near-1 h² or REML SE ≥ 0.3, and aggregate edu/birth-year
summaries per stratum.

### Cross-cell h² heterogeneity test

`gcta --reml --mgrm` cannot run for disjoint-sample stratum GRMs (it
segfaults; see Plan_deviations §7b for the full diagnosis). The test
implemented here is **Cochran's Q** on the 4 per-cell h² estimates:

```
Q = Σ_k (h²_k − h²_pooled)² / SE(h²_k)²
h²_pooled = Σ_k w_k h²_k / Σ_k w_k,   w_k = 1 / SE(h²_k)²
Q ~ χ²_{K-1} under H₀: all h²_k equal.
```

Plus pairwise Z-tests `Z_ij = (h²_i − h²_j) / √(SE_i² + SE_j²)` and the
I² heterogeneity statistic. Output: `output/final/heterogeneity.tsv`.

### Power calc (not reported in headline; diagnostic only)

The Visscher (2014, PLOS Genet 10(4): e1004269) approximation
`SE(h²) ≈ √(2 / (N² · Var(off-diag-GRM)))` is implemented in
`scripts/08_power_calc.R` and evaluates a per-cell SE/power across the
`H2_GRID` from the empirical `Var(off-diag)` of each `*_unrel` GRM.

**These values are not reported in the headline.** Analytical power
estimates based on empirical GRM variance are invalid for the
0.20-pruned cell-level analysis because the retained sample contains
substantial relatedness and violates the unrelated-sample assumption
underlying the Visscher approximation. Empirically, the Visscher SE
fell well below the REML SE in the two cells that retained the most
relatedness at cutoff 0.20, implying near-certain power that the REML
standard errors contradict. **Precision is therefore summarized using the
GCTA REML SE** from each `.hsq` file. The `power.tsv` file is generated
Tardis-side by `scripts/08_power_calc.R`, copied to
`output/final/power.tsv` by step 09, transferred to the Mac shared
volume by `transfer_results_from_tardis.sh`, and is **not committed to
the git repo** (consistent with the repo-wide policy that data
products live on the shared volume, not in version control; see
`.gitignore`).

**Scope.** `power.tsv` is computed for the pooled-GRM-sensitivity track
only (`analysis = main_4cell`). It is not computed for
the block-diagonal primary track because the empirical Var(off-diag)
on a block-diagonal GRM is dominated by zero cross-cohort entries by
construction, which would make the Visscher SE artificially large and
the resulting power numbers uninformative. They are not computed for
the per-cohort sensitivity or the R1/G1/RG1 secondary tracks either;
precision in those tracks is summarised by the REML SE per stratum in
the corresponding TSVs.

---

## Methodological diagnostics

### Cross-cohort GRM cohort-fingerprint check

The pooled GRM combines 5 cohorts that differ in imputation panel
(HRC vs 1000G) and array platform. In principle, pooled-allele-frequency
standardization in GRM construction can turn real cross-cohort
allele-frequency differences into spurious off-diagonal relatedness —
"cohort fingerprint" — which inflates V_G and biases h² upward.

A diagnostic partitions each unrel GRM's off-diagonal entries into
WITHIN-cohort vs ACROSS-cohort pairs and compares their means. A
systematic excess of WITHIN over ACROSS would indicate residual
fingerprint that the `--grm-cutoff 0.20` filter didn't absorb.

The diagnostic, run on the `--grm-cutoff 0.05` unrelated GRMs, did **not**
detect a large mean shift between within-cohort and across-cohort GRM
entries (cohort fingerprint is even less of a concern at the higher cutoff
the pooled-GRM sensitivity uses). It does not rule out residual cohort
structure in the GRM distribution, tail behaviour, or eigenstructure, nor
does it assess cohort predictability from GRM components or whether
phenotype differences across cohorts load onto genetic structure. The
pooled-GRM sensitivity is therefore interpreted alongside the primary
block-diagonal track (steps `03d/04d/05d/06d`, with cohort fixed effects),
which by construction sets cross-cohort GRM entries to zero.

The *overall* mean off-diagonal varies across strata and is highest in
East × post-1990. This is consistent with residual regional relatedness
or population structure in the SHIP-dominated Pomeranian sample. It is
not a technical cohort fingerprint, but it is also not necessarily
innocuous: regional genetic structure can align with environmental and
educational structure, which complicates interpretation of
within-stratum h². The elevated off-diagonals also increase the
empirical Var(off-diag-GRM), which is what made the Visscher (2014)
power formula spuriously confident in that cell (see *Methods → Power
calc*).

### Cohort × cell retention at `--grm-cutoff`

The cohort-fingerprint diagnostic above tests the *mean* of off-diagonals.
A separate concern is whether the *tail* (pairs > cutoff) is balanced
across cohorts inside each cell — i.e., is `--grm-cutoff` dropping
minority-cohort individuals more aggressively than dominant-cohort ones?
A composition diagnostic checks this by comparing the analytic-sample
composition to the post-cutoff GRM composition, broken down by
cohort × cell.

At `GRM_CUTOFF=0.05` on the pooled cross-cohort GRM, the minority cohort
in each cell was eliminated at >99%, silently restricting the h² estimates
to a non-representative subsample; that is the reason the pooled-GRM
sensitivity runs at 0.20, where retention is uniform (~80–100%) across
cohorts within each cell. See `Plan_deviations.md` §7e for the per-cell
retention tables.

---

## Output files

```
output/
├── geno/                                # per-cohort INFO≥.90 bfiles + pooled
│   ├── BASEII_harmonized_r2_0.9.{bed,bim,fam}
│   ├── ... (one per cohort) ...
│   └── pooled.{bed,bim,fam}
├── grm/                                 # per-stratum + per-cohort + block-comp GRMs
│   ├── grm_<cell>.grm.{bin,id,N.bin}                      # pooled-GRM sensitivity (steps 03/04)
│   ├── grm_<cell>_unrel.grm.{bin,id,N.bin}                # post-cutoff 0.20 (step 04)
│   ├── grm_cohort_<COHORT>.grm.{bin,id,N.bin}             # per-cohort, pre-cutoff (step 03c)
│   ├── grm_cohort_<COHORT>_unrel.grm.{bin,id,N.bin}       # per-cohort, post-cutoff 0.05 (step 04c)
│   ├── grm_blockcomp_<COHORT>_<STRATUM>.grm.{bin,id,N.bin} # cohort × stratum blocks (step 03d for cells; 03e for R1/G1/RG1 strata)
│   └── grm_blockdiag_<STRATUM>.grm.{bin,id,N.bin}         # combined block-diagonal GRMs (step 04d); STRATUM ∈ {<cell>, region_*, gender_*, region_gender_*}
├── reml/                                # .hsq files + Q/Z + power
│   ├── main/standalone_*.hsq            # pooled-GRM sensitivity: 8 fits (4 cells × 2 variants)
│   ├── per_cohort/*.hsq                 # per-cohort sensitivity (step 05c)
│   ├── blockdiag/*.hsq                  # cell-level primary (step 05d, default)
│   ├── blockdiag_region_main/*.hsq      # R1 secondary track (step 05d, --step xtra)
│   ├── blockdiag_gender_main/*.hsq      # G1 secondary track
│   ├── blockdiag_region_gender/*.hsq    # RG1 secondary track
│   ├── heterogeneity.tsv                # step 07 output (pooled-GRM sensitivity)
│   ├── lrt.tsv -> heterogeneity.tsv     # compatibility symlink
│   └── power.tsv                        # step 08 output (diagnostic only)
└── final/                               # ← deliverables, rsync these back to the Mac
    ├── h2_snp.tsv                       # pooled-GRM sensitivity table
    ├── heterogeneity.tsv                # pooled-GRM Q + pairwise Z
    ├── power.tsv                        # Visscher diagnostic (not headline)
    ├── manifest.tsv
    ├── dashboard_row.md
    ├── h2_blockdiag.{tsv,md}            # cell-level primary
    ├── h2_blockdiag_heterogeneity.tsv
    ├── h2_per_cohort.{tsv,md}           # per-cohort sensitivity
    ├── h2_region_main.{tsv,md}          # R1 secondary track
    ├── h2_region_main_heterogeneity.tsv
    ├── h2_gender_main.{tsv,md}          # G1 secondary track
    ├── h2_gender_main_heterogeneity.tsv
    ├── h2_region_gender.{tsv,md}        # RG1 secondary track
    └── h2_region_gender_heterogeneity.tsv
```

Plus, after Mac-side `summarize_results.R`:

```
${DATA_ROOT}/04_GxE_Analysis/02_SNP_Heritability/output/final/
├── summary_report.md                    # human-readable, embeds key numbers
├── h2_forest_plots.pdf                # combined multi-page PDF: one forest plot per analysis
│                                      #   page 1: primary cell-level (block-diagonal + cohort FE)
│                                      #   page 2: pooled-GRM sensitivity (cell-level, cutoff 0.20)
│                                      #   page 3: per-cohort sensitivity
│                                      #   page 4: R1 Region main
│                                      #   page 5: G1 recorded gender/sex main
│                                      #   page 6: RG1 Region × gender/sex
├── h2_forest_plot.png                 # PNG: pooled-GRM sensitivity (for inline embed)
├── h2_blockdiag_forest_plot.png       # PNG: primary cell-level (for inline embed)
├── h2_per_cohort_forest_plot.png      # PNG: per-cohort sensitivity
├── h2_region_main_forest_plot.png     # PNG: R1
├── h2_gender_main_forest_plot.png     # PNG: G1
└── h2_region_gender_forest_plot.png   # PNG: RG1
```

---

## Limitations

- **The primary block-diagonal model adjusts cohort means but does not model cohort-specific residual variances or cohort-specific V_G.** `--covar gender_cohort.covar` removes between-cohort mean differences in years of education as fixed effects, but the single-V_G / single-V_e GREML specification still imposes shared variance components across cohort blocks. If cohorts differ substantially in residual variance or in within-cohort h², the primary estimate is a compromise across blocks; the diagnostic in `Plan_deviations.md` §7g indicates that for the more informative pre-1990 cells in this sample, the omitted-cohort-mean problem was the dominant misspecification, but heterogeneous within-cohort residuals are not formally modelled. The within-cohort Var(education) max/min ratio inside each primary cell is ≤ 1.44 (East × pre-1990 1.44; East × post-1990 1.41; West × pre-1990 1.17; West × post-1990 1.02) — i.e. the within-cell variance ratios do not indicate severe cohort-level heterogeneity in education dispersion. The single-V_e assumption remains a modelling simplification, but it is not strongly contradicted by these diagnostics.
- **Post-1990 cells remain imprecise.** Across specifications, the post-1990 cells were consistently imprecise, with confidence intervals covering most of the admissible range. Boundary or near-boundary point estimates in those cells are boundary fits from imprecise models and are not informative for cell-specific inference. These cells contribute as data points to Cochran's Q but do not support cell-specific inference.
- **East × pre-1990 remains SHIP-heavy.** In the pooled-GRM sensitivity at cutoff 0.20, the iterative `--grm-cutoff` removal silently excludes most SOEP-East and TwinLife-East participants, leaving predominantly SHIP plus BASE-II East Berlin. In the primary block-diagonal track, the SHIP block dominates the East × pre-1990 cell's information by sample size as well. Either way, the East × pre-1990 estimate is not a balanced estimate over all East-pre Germans; it is closest to a SHIP-Pomerania-weighted estimate.
- **Pooled-GRM sensitivity at cutoff 0.20** retains substantial relatedness (it is not an unrelated-sample GREML) and carries potential residual cross-cohort allele-frequency structure even though the within–across mean diagnostic does not detect a large fingerprint. Reported as a sensitivity, not the primary.
- **Per-cell SNP set may differ across cells in the pooled-GRM sensitivity** because `gcta --maf 0.01` is applied within each cell's `--keep` subset. The primary block-diagonal track inherits the cohort-level SNP set + within-cohort standardization from the per-cohort GRMs and is more uniform.
- **SOEP chr1 INFO file was never delivered by the data provider.** `03a_filter_imputation_quality.sh` falls back to including all chr1 variants from the SOEP harmonized `.bim` at the mildQC R² > 0.1 threshold. Reflected in *Methods → Variant set*, *TL;DR → SNP set, by track*, and in `manifest.tsv` as `soep_chr1_info=fallback_mildQC_R2_0.1`. A sensitivity excluding SOEP chr1 fallback variants is not implemented.
