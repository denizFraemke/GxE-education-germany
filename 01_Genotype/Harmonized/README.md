# PGI Computation Pipeline — Cross-Cohort Harmonized

Polygenic Index (PGI) computation for the project
**"Genomics of Educational Attainment Across 80 Years of Social and Political Transformation in Germany"**.

Five German cohorts, four traits (EA4 / Cog / NonCog / Height), SBayesR weights, a fully harmonized cross-cohort SNP set, and final TSV tables ready for downstream GxE analysis.

- 5 cohorts, N = 18,898 total, scored on the intersection of 6,878,689 shared SNPs (R² ≥ 0.3 across all cohorts).
- All scripts read their settings from `config.sh`; the only command-line arguments are `run_pipeline.sh --from 00` / `--step NN`.

---

## Reported PGI validation

| Trait  | Result |
|--------|--------|
| EA4    | Pearson r with education: BASEII 0.24 (R² = 0.059), SHIP pooled 0.29 (R² = 0.081), SOEP 0.25 (R² = 0.063), TwinLife 0.34 (R² = 0.118) |
| Cog    | Scored in all five cohorts; no leave-one-out weight needed (none of the cohorts is in the Cog discovery GWAS) |
| NonCog | Scored in all five cohorts with the published NCog weights |
| Height | Negative control, r ≈ 0.35 |

Reported correlations are Pearson r between `PGI_EA4_global_z` and the standardized education measure in the phenotype file; sample sizes: BASEII n = 1,410; SHIP pooled n = 7,906; SOEP n = 2,425; TwinLife n = 3,331.

---

## How to run this

The pipeline runs in three explicit phases: local setup on a Mac with the
institution SMB shares mounted, the SLURM compute on the cluster, and pulling
the deliverables back to the downstream analysis project. Phases 1 and 3 use
the Mac as a transit machine; phase 2 does all the work on the cluster.

### Phase 1 — Local setup (Mac, SMB volumes mounted)

**Step 0 — one-time per-user setup.** Every cluster/SSH-touching script
reads its identity from `user_config.sh` (gitignored). The first
`transfer_*.sh` you run walks you through writing it interactively
(username, project root, bin dir, LD-reference dir, and the local
analysis-project path).

Then verify SSH-key auth to the cluster (the transfer scripts run rsync
non-interactively, so password prompts will hang them):

```bash
source user_config.sh    # loads $TARDIS_USER, $TARDIS_HOST, ...
ssh "$TARDIS_USER@$TARDIS_HOST" echo ok   # must NOT prompt for a password
# If it does:  ssh-copy-id "$TARDIS_USER@$TARDIS_HOST"
```

**Then stage data + code.** Mount the five project shares
(Finder → Go → *Connect to Server*,
`smb://<institution-share>/<volume>`); they appear under
`/Volumes/<volume>/`. From inside `01_Genotype/Harmonized/`:

```bash
# 1. Push BASEII / SHIP-0 / SHIP-Td / SOEP .info.gz files to the cluster under
#    ${TARDIS_PROJECT_ROOT}/data/info_files/. (TwinLife has no .info.gz —
#    its R² TSV rides along in step 2 below.)
bash scripts/transfer_info_files.sh
#    Alternative on Posit Workbench / ARC (no Mac SMB, much faster):
#    bash scripts/extract_and_transfer_info_workbench.sh

# 2. Stage cohort genotype bfiles + the TwinLife imputation-R² TSV to the
#    cluster. macOS:
bash transfer_geno_to_tardis.sh
#    On Posit Workbench / ARC (Linux, uses biomount instead of SMB):
#    bash transfer_geno_to_tardis_workbench.sh
#    Both scripts: bfiles → ${TARDIS_PROJECT_ROOT}/data/geno/,
#                  TwinLife_imputation_r2.tsv.gz →
#                    ${TARDIS_PROJECT_ROOT}/data/info_files/TWINLIFE/
#    Resume via *.completed state files and retry on SMB drop-outs.

# 3. Stage the formatted GWAS summary statistics (per-chr SBayesR input) to the
#    cluster under ${TARDIS_PROJECT_ROOT}/data/formatted_sumstats/. They live on
#    the data share at Projects/04_data_analysis/013-PGS/data/formatted/
#    (one subdir per trait; names match config.sh TRAIT_SUMSTATS_SUBDIR).
rsync -av \
  "${SHARE_ROOT}/Projects/04_data_analysis/013-PGS/data/formatted/"{EA4_excl_BASEII,EA4_excl_SHIP,Cog_Malanchini,NCog_Malanchini,Height_GIANT} \
  "$TARDIS_USER@$TARDIS_HOST:$TARDIS_PROJECT_ROOT/data/formatted_sumstats/"
#    (Raw Malanchini Cog/NCog GWAS sit alongside in ../sumstats/.)

# 4. Push the pipeline code (and your user_config.sh) to the cluster.
#    Excludes output/, logs/, data/. user_config.sh IS synced — that's
#    how config.sh on the cluster picks up your bin / LD-ref paths.
bash transfer_pipeline_to_tardis.sh
```

`${SHARE_ROOT}` is the mount point of the institutional data share (default
`/path/to/share_root`). It is distinct from `${DATA_ROOT}`, which the merge and
analysis stages use for the project data root inside that share.

The TwinLife R² TSV is a pre-extracted aggregate file deposited on the
TwinLife share at
`${TWINLIFE_SHARE}/private/data/Gendata/rawgendata/TwinLife_imputation_r2.tsv.gz`
(everyone with TwinLife data access can read it). The transfer scripts
above pick it up automatically from there.

If the SOEP-G genotypes were delivered as per-chromosome GPG-encrypted
files (the typical case), decrypt and merge them once on the cluster after
the transfer:

```bash
ssh "$TARDIS_USER@$TARDIS_HOST"
cd "$TARDIS_PROJECT_ROOT"
bash scripts/setup_soep_geno.sh --passphrase-file ~/soep_passphrase.txt
```

See [Cohort notes](#cohort-notes) for the rest of the per-cohort quirks.

### Phase 2 — On the cluster (SLURM)

The pipeline splits into two halves: SBayesR (a SLURM array, hours of
wallclock) and steps 02–07b (login-node or a single SLURM job, ~1–2 h
total). Run them as two steps (useful when SBayesR fails mid-array, or
when you want SLURM to chain the second half via a dependency):

```bash
# (a) Submit SBayesR; run_pipeline.sh stops automatically after step 01.
bash run_pipeline.sh --from 00 --to 01

# (b) Wait until the sbr_* jobs leave the queue.
squeue -u $USER | grep sbr_

# (c) Continue. As a SLURM job (recommended for long sessions):
sbatch submit_post_sbayesr.slurm
#     ...or interactively on the login node:
bash run_pipeline.sh --from 02
```

Step 04c (ancestry keep-lists) is run separately — see the step table below.

### Phase 3 — Pull results back

Two directories on the cluster hold everything the downstream analysis
project consumes. The convenience wrapper pulls both (PGI tables, PC
eigenvalue tables for the scree, summaries, diagnostics) to the analysis
project's `data/final/` on the SMB share:

```bash
bash 01_Genotype/Harmonized/transfer_final_from_tardis.sh
```

The equivalent manual `rsync` (set `ANALYSIS_ROOT` in `user_config.sh`, or
override the env var for this one command):

```bash
# From the Mac, with SMB volumes mounted:
source 01_Genotype/Harmonized/user_config.sh   # loads $TARDIS_SSH, $TARDIS_PROJECT_ROOT, $ANALYSIS_ROOT

rsync -av "${TARDIS_SSH}:${TARDIS_PROJECT_ROOT}/output/final/"       "${ANALYSIS_ROOT}/data/final/"
rsync -av "${TARDIS_SSH}:${TARDIS_PROJECT_ROOT}/output/diagnostics/" "${ANALYSIS_ROOT}/data/final/diagnostics/"
```

What lands locally:

- `output/final/{COHORT}_PGI_PCs.tsv` — five per-cohort deliverable tables (the analysis project consumes these).
- `output/final/{COHORT}_PC_eigenval.tsv` — five per-cohort PCA eigenvalue tables (`cohort, PC, eigenval`), emitted by step 04. `03_Merge/validate.R` reads these to draw the ancestry-PC scree (variance explained per PC). They live in `output/final/` alongside the `_PGI_PCs.tsv` files, so the same `rsync` pulls them — keep them in the transfer.
- `output/final/PGI_summary_statistics.tsv` — summary table.
- `output/diagnostics/global_pgi_statistics.tsv` and the rest of the diagnostics tree — QC TSVs.

The full `output/` tree is large and regenerable; do not mirror it.

---

## Pipeline DAG

Every numbered script that participates in the production DAG, in run
order. Wall-times are rough order-of-magnitude estimates measured on the
cluster with the 5-cohort inputs.

| Step | Script | Phase | Input | Output | Typical wall-time |
|------|--------|-------|-------|--------|-------------------|
| 00   | `scripts/00_check_inputs.sh`              | 2 (cluster) | `data/geno/`, `data/formatted_sumstats/`, `data/info_files/`, `${LD_REF_DIR}` | (log only) | < 1 min |
| 01   | `scripts/01_submit_all_sbayesr.sh` → `scripts/01_run_sbayesr.slurm` (SLURM array 1–22) | 2 (cluster SLURM) | `data/formatted_sumstats/<TRAIT>/*.ma`, LD matrices | `output/sbayesr/<TAG>/<TAG>_chr{1..22}_sbayesR.snpRes` | ~1–6 h wallclock per (trait × chr) array; full SBayesR pass ≈ 6 h end-to-end |
| 02   | `scripts/02_concat_weights.sh`            | 2 (cluster) | per-chr `.snpRes` from step 01 | `output/weights/<TAG>_sbayesR_chrpos.txt` | ~5 min |
| 03a  | `scripts/03a_filter_imputation_quality.sh` | 2 (cluster) | `data/info_files/<COHORT>/*.info.gz` (+ TwinLife R² TSV) | `output/geno/<COHORT>_r2_pass_{0.3,0.8}.txt` | ~10 min |
| 03   | `scripts/03_harmonize_geno.sh`            | 2 (cluster) | `data/geno/<COHORT>.{bed,bim,fam}` + R² pass lists | `output/geno/<COHORT>_harmonized_shared.{bed,bim,fam}` (6,878,689 shared SNPs) | ~30 min |
| 03b  | `scripts/03b_fix_strand.sh` *(off by default)* | 2 (cluster) | harmonized bfiles | flipped bfiles | n/a (disabled — `FIX_STRAND=no`) |
| 03c  | `scripts/03c_align_to_ukb.sh` *(off by default)* | 2 (cluster) | harmonized bfiles + UKB reference | UKB-aligned bfiles | n/a (disabled — `UKB_ALIGN=no`) |
| 04   | `scripts/04_compute_pcs.sh`               | 2 (cluster) | harmonized bfiles | `output/pcs/<COHORT>_pcs.eigenvec` | ~30 min per cohort |
| 04b  | `scripts/04b_compute_cross_cohort_pcs.sh` | 2 (cluster) | merged harmonized bfiles | `output/pcs/cross_cohort_pcs.eigenvec` | ~45 min |
| 04c  | `scripts/04c_ancestry_1kg.sh` → `scripts/04c_classify_eur.R` *(run separately, not driven by `run_pipeline.sh`)* | 2 (cluster) | harmonized bfiles + 1000G phase-3 b37 reference VCFs | `output/ancestry/keep_lists/EUR_keep_IDs_<KEY>.tsv` — the European-ancestry keep-lists `03_Merge` consumes via `filter_to_eur()` | dominated by the one-time 1000G reference build |
| 05   | `scripts/05_score_pgi.sh`                 | 2 (cluster) | harmonized bfiles + weights | `output/scores/<COHORT>_<trait>.sscore` | ~15 min |
| 06b  | `scripts/06b_r2_sensitivity.sh`           | 2 (cluster) | strict-R² shared SNP set + harmonized bfiles | `output/scores/<COHORT>_<trait>_r2strict.sscore` | ~15 min |
| 07   | `scripts/07_assemble_output.R`            | 2 (cluster) | step 05 / 06b sscores + step 04 / 04b eigenvecs | `output/final/<COHORT>_PGI_PCs.tsv` + global stats | ~5 min |
| 07b  | `scripts/07b_compare_harmonization.R`     | 2 (cluster) | `output/final/*` | `output/diagnostics/harmonization_*.tsv` | ~5 min |

Step 04c has its own method write-up in
[`scripts/README_04c_ancestry.md`](scripts/README_04c_ancestry.md), including
the public reference download and the classification parameters.

Numbered scripts that are **not** in the production DAG (kept in `scripts/` for one-off use; see their own headers):

- `scripts/11_ea4_harmonization_loss.sh` — cluster-side diagnostic that
  attributes EA4 SNP loss across the 5-cohort intersection.
- `scripts/08_plot_pgi_qc.R`, `scripts/09_compare_education_ea4_same_sample.R`
  — local Mac-side QC plots and cross-checks against the downstream analysis
  project's merged RDS.

---

## Cohort notes

Quirks that affect how each cohort moves through the pipeline.

- **BASE-II** (Affymetrix 6.0 / HRC v1.1, n = 2,351). Standard HRC
  pipeline.
- **SHIP-0** (Affymetrix 6.0 / HRC v1.1, n = 4,069). Standard HRC
  pipeline.
- **SHIP-Td** (Illumina Omni 2.5 + GSA / HRC v1.1, n = 4,120). A
  post-imputation merge of two genotyping batches imputed separately.
  Step 03a takes `min(R²_batch1, R²_batch2)` per variant — a SNP must
  be well-imputed in **both** batches to pass the threshold.
- **SOEP-G** (Illumina GSA / HRC v1.1, n = 2,497). Delivered as 66
  GPG-encrypted per-chromosome VCFs (the SOEP data-sharing agreement
  rules out unencrypted transit). On the cluster, run
  `scripts/setup_soep_geno.sh --passphrase-file ...` once after the
  transfer to decrypt + merge into `data/geno/SOEP-G.b37.ext_chrall.{bed,bim,fam}`.
  The SOEP `chr1.info.gz` is not part of the delivery: the SOEP info-file
  directory on the source share
  (`${SOEP_SHARE}/private/data/001_SOEP_rawdata/004_SOEP_SNP_level/soep-g_ext_share/info_files/`)
  contains `chr2..chr22` only. Both transfer scripts skip the missing chr1
  cleanly — `scripts/transfer_info_files.sh` uses `chr*.info.gz` globbing
  (and probes chr2 for the pre-flight existence check), and
  `scripts/extract_and_transfer_info_workbench.sh` iterates `chr1..22` with
  per-file `[[ -f ]]` guards. `scripts/03a_filter_imputation_quality.sh`
  detects the missing chr1 info file at runtime and falls back to including
  **all** chr1 variants from the SOEP harmonized `.bim`, which inherits the
  provider's mildQC R² > 0.1 filter (see the `mildQC` filename suffix). Net
  effect: SOEP chr1 variants on any downstream pass list are at R² > 0.1,
  not the documented threshold. They survive the cross-cohort `comm -12`
  intersection only if the other four cohorts pass them at the documented
  threshold, so the worst case is a mild precision cost on chr1 rather than
  bias.
- **TwinLife** (Illumina ILL_TWL / **1000G Phase 3 v5**, n = 5,861).
  The odd cohort out: 1000G panel instead of HRC, and no `.info` files
  in the raw delivery — R² is embedded in the dosage VCF INFO field.
  The extracted per-SNP R² values are deposited as
  `TwinLife_imputation_r2.tsv.gz` on the TwinLife share alongside
  the dosage VCFs and ride along with the genotype transfer
  (`transfer_geno_to_tardis.sh`) → `data/info_files/TWINLIFE/` on the cluster.
  TwinLife is also the SNP-count bottleneck (~10.8 M imputed variants
  vs ~39 M for the HRC cohorts); the 6.88 M shared intersection is
  driven by its catalogue.

All five cohorts are scored with `EA4_excl_SHIP` for EA4 and with the
published Malanchini `NCog` for NonCog; `Cog` and `Height` are the same
weights everywhere. See [Leave-One-Out weight mapping](#leave-one-out-weight-mapping-full-matrix).

---

## Data Inventory

### Cohorts

| Cohort   | N     | Genotyping array              | Imputation panel   | Raw n_variants | After R² ≥ 0.3 | Shared set |
|----------|-------|-------------------------------|--------------------|----------------|----------------|------------|
| BASE-II  | 2,351 | Affymetrix 6.0                | HRC v1.1           | ~39.1 M        | ~12 M          | 6,878,689  |
| SHIP-0   | 4,069 | Affymetrix 6.0                | HRC v1.1           | ~39.1 M        | ~12 M          | 6,878,689  |
| SHIP-Td  | 4,120 | Illumina Omni 2.5 + GSA (merge) | HRC v1.1         | ~37.2 M        | ~12 M          | 6,878,689  |
| SOEP-G   | 2,497 | Illumina GSA (GSA2019)        | HRC v1.1           | 23,185,408*    | ~12 M          | 6,878,689  |
| TwinLife | 5,861 | Illumina (ILL_TWL)            | **1000G Phase 3 v5** | ~10.8 M      | ~7.6 M         | 6,878,689  |
| **Total**| **18,898** | — | — | — | — | 6,878,689  |

**TwinLife is the SNP-count bottleneck.** It uses 1000G Phase 3 instead of HRC v1.1, producing ~10.8 M imputed variants vs ~39 M for the HRC cohorts. The shared intersection (~7.6 M pre-R² / 6.88 M post-R²) is therefore driven by TwinLife's smaller variant catalog.

\* SOEP is delivered already filtered by the provider (filename carries `mildQC`). The 23.2 M is the merged-bim count on the cluster (`wc -l SOEP-G.b37.ext_chrall.bim`); the per-chr Minimac info-file sum is 20,461,505, so the bim includes ~2.7 M multiallelic splits.

### Genotyping and imputation provenance

| Detail | BASE-II | SHIP-0 | SHIP-Td (merged) | SOEP-G | TwinLife |
|---|---|---|---|---|---|
| **Genotyping array** | Affymetrix Genome-Wide Human SNP Array 6.0 | Affymetrix Genome-Wide Human SNP Array 6.0 | Batch 1: Illumina HumanOmni 2.5; Batch 2: Illumina GSA (GSAMD-24v1) | Illumina GSA (GSA2019) | Illumina ILL_TWL |
| **Calling algorithm** | not documented in the delivery | Birdseed2 | Batch 1: GenomeStudio GenCall; Batch 2: GenomeStudio 2.0 GenCall | not documented in the delivery | not documented in the delivery |
| **Original build** | b36 (NCBI Build 36), lifted to b37 | b36 (hg18), forward-strand alleles | Batch 1: b36 (hg18); Batch 2: GRCh37 (hg19) | GRCh37 (b37) | GRCh37 (b37), `human_g1k_v37.fasta` reference |
| **Build in the pipeline** | GRCh37/b37 | GRCh37/b37 (lifted during imputation) | GRCh37/b37 (lifted/merged during imputation) | GRCh37/b37 | GRCh37/b37 |
| **Imputation reference panel** | HRC | **HRC v1.1** | **HRC v1.1** (both batches) | **HRC v1.1** | **1000 Genomes Phase 3 v5** |
| **Imputation software** | Minimac3 | Michigan Imputation Server (Eagle v2 + Minimac3) | Michigan Imputation Server (Eagle v2 + Minimac3) | Michigan Imputation Server (Minimac) | Minimac4 v1.0.2 |
| **Imputation date** | 2019-09-26 | 2016-07-29 | Batch 1: 2016-07-29; Batch 2: 2017-11-17; merged: 2018-10-25 | not documented (VCFs are GPG-encrypted) | 2022-09-24 |
| **SNPs pre-imputation** | not documented in the delivery | 823,635 | Batch 1: 1,803,558; Batch 2: 488,725 | not documented in the delivery | 745,980 |
| **SNPs post-imputation** | ~39,131,578 | ~40,356,094 (39,127,678 in the final bed) | Batch 1: ~40,360,311; Batch 2: ~40,359,613 (37,183,531 in the merged bed) | not documented (encrypted source) | ~10,832,804 |
| **Pre-imputation QC by the provider** | not documented in the delivery | pHWE > 0.0001, call rate > 0.95, non-monomorphic; b36→b37 liftover; duplicates and inconsistent-reference variants removed | Batch 1: pHWE > 0.0001, call rate > 0.95, non-monomorphic; Batch 2: additionally MAC ≥ 10, MAF ≥ 1%; both lifted to b37 | "mildQC" (filename; details not documented) | not documented in the delivery |
| **Post-imputation R²/INFO filter by the provider** | none | none | none (both batches) | not documented | not documented |
| **Strand orientation** | forward | forward | forward | forward | forward |
| **R² info files delivered** | yes (`LIFEBRAIN_BASEII_2015_AFFY.info.gz`) | yes (`SHIP-0_R4a.chr*.info.gz`) | yes (`SHIP-Td.chr*.info.gz` + `SHIP-Td_B2.chr*.info.gz`) | yes (`SOEP-G.b37.mildQC.hrc1-1_imp.chr*.info.gz`, chr2..22) | no (R² is embedded in the dosage VCF INFO field) |
| **Dosage data delivered** | yes (`.dose.vcf.gz`) | yes (`.dose.vcf.gz`) | yes (merged `.vcf.gz` with DS field) | yes (`.dose.vcf.gz`, GPG-encrypted) | yes (`*.imputed.dose.vcf.gz`) |
| **Pipeline input format** | bed/bim/fam hard-call, converted from VCF | bed/bim/fam hard-call, converted from VCF | bed/bim/fam hard-call, converted from the merged VCF | bed/bim/fam hard-call, decrypted + converted | bed/bim/fam hard-call, converted from VCF |
| **Data custodian** | LIFEBRAIN consortium | University of Greifswald (SHIP) | University of Greifswald (SHIP), merged there and converted by the data provider | SOEP-G (DIW Berlin), GPG-encrypted delivery | TwinLife study |

Three array families are represented — Affymetrix 6.0 (~820 K genotyped SNPs),
Illumina Omni 2.5 (~1.8 M) and Illumina GSA (~490–750 K) — and two imputation
panels (HRC v1.1 for four cohorts, 1000G Phase 3 v5 for TwinLife). The merged
SHIP-Td file carries fewer variants (37.2 M) than either batch alone (~40.4 M)
because the merge intersects them. No provider applied a post-imputation R²
filter, which is why step 03a applies one.

### Reference data

- **LD reference** for SBayesR: UKB 50k European, shrunk sparse matrices (Zeng et al.), ~2.8 M SNPs, MAF > 0.01, per-chromosome.
- **GWAS summary statistics**:
  - EA4 (Okbay et al., 2022): two leave-one-out versions are prepared, `EA4_excl_BASEII` and `EA4_excl_SHIP`. Scoring uses `EA4_excl_SHIP` for all five cohorts (see [Weight mapping](#weight-mapping-single-loo-weight-file)).
  - Cog / NCog (Malanchini et al., 2024): the published weights are used for all five cohorts. See `Plan_deviations.md` §22 for why the study-excluded NonCog is not used.
  - Height (Yengo et al., 2022 / GIANT): negative control. **The pre-formatted Height copy on the data share is malformed** — it keeps GIANT's leading `CHR` column (9 cols instead of 8) and leaves in SNPs with missing effect estimates (NA). GCTB ≥ 2.5 reads those columns shifted, matching 0 SNPs and segfaulting on the first NA row. Regenerate the correct 8-column `.ma` from the **raw** GIANT file (also on the share: `Projects/03_data/009_SUMSTATS/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz`) with [`scripts/format_height_sumstats.sh`](scripts/format_height_sumstats.sh) — set `GIANT_HEIGHT_RAW` to that path. Step 00 verifies the 8-column layout and points here if it isn't.
- **Where the formatted sumstats are staged from**: the pre-formatted per-chromosome SBayesR files live on the data share at `Projects/04_data_analysis/013-PGS/data/formatted/` (subdir names match the trait names above). `scripts/setup_tardis_data.sh` copies/symlinks them into `data/formatted_sumstats/` on the cluster. **Exception — `Height_GIANT/`: do _not_ stage it from there (that copy is malformed, see above); regenerate it from raw GIANT with `scripts/format_height_sumstats.sh`.** The repo therefore stays reproducible from share inputs alone: the other traits are copied pre-formatted, Height is rebuilt from the raw GIANT file.

### Tool versions

- PLINK v2.0.0-a.7LM (7 Mar 2026)
- PLINK v1.90b7.1 (18 Oct 2023)
- GCTB 2.5.2 (SBayesR)

---

## EA4 — the primary outcome

EA4 is the educational attainment PGI and the scientifically central score in this project. Everything about how it is computed deserves to be explicit.

### Weight mapping (single LOO weight file)

All five cohorts are scored with the `EA4_excl_SHIP` leave-one-out weight file:

| Cohort       | EA4 weight file       | Rationale |
|--------------|-----------------------|-----------|
| BASE-II      | `EA4_excl_SHIP`       | BASE-II contributes ~0.07 % of the EA4 discovery sample, so the in-sample bias from scoring it against `EA4_excl_SHIP` is statistically undetectable (r² inflation ≈ 1.0007) |
| SHIP-0       | `EA4_excl_SHIP`       | SHIP is in discovery — excluded |
| SHIP-Td      | `EA4_excl_SHIP`       | SHIP is in discovery — excluded |
| SOEP-G       | `EA4_excl_SHIP`       | Not in discovery |
| TwinLife     | `EA4_excl_SHIP`       | Not in discovery |

A single weight file across all five cohorts puts every PGI on the same scale
without an analytical or empirical calibration step, so cross-cohort raw-mean
differences reflect population variation rather than LOO methodology.

The mapping is enforced centrally in `config.sh` by `get_ea4_trait_for_cohort()` and consumed by `scripts/05_score_pgi.sh`. This helper is the **single source of truth**; the same applies to `get_noncog_trait_for_cohort()` for NonCog.

### SBayesR parameters (EA4)

- Gamma mixture: `0, 0.01, 0.1, 1`
- Pi: `0.95, 0.02, 0.02, 0.01`
- Chain length 10,000, burn-in 2,000
- MHC excluded by default
- LD reference: UKB 50k shrunk sparse per chromosome

The SBayesR output is converted to a three-column chr:pos / A1 / BETA scoring file by `02_concat_weights.sh` and consumed by PLINK2 `--score`.

### Final validation

After running steps 00 → 07:

| Cohort       | n    | Pearson r(PGI_EA4_global_z, education) |
|--------------|------|----------------------------------------|
| BASE-II      | 1,410 | 0.242 |
| SHIP (pooled)| 7,906 | 0.285 |
| SOEP-G       | 2,425 | 0.250 |
| TwinLife     | 3,331 | 0.343 |

Cog / NonCog / Height diagnostics are in `output/final/PGI_summary_statistics.tsv`.

---

## Cross-cohort harmonization

Six layers of harmonization ensure that PGIs are comparable across cohorts despite different genotyping platforms and imputation panels.

1. **Variant-ID normalization** (`03_harmonize_geno.sh`): IDs set to `chr@:#` format via PLINK2 `--set-all-var-ids`; duplicate positions resolved with `--rm-dup force-first`.
2. **Imputation-quality filter** (`03a_filter_imputation_quality.sh`): per-cohort R² pass lists at threshold 0.3 (default) and 0.8 (sensitivity). Info files are in `data/info_files/` (Minimac format; TwinLife uses a custom tsv).
3. **Shared-SNP intersection**: after variant-ID normalization + R² filter, all five cohorts are restricted to the intersection of their variant sets → **6,878,689 shared SNPs**, byte-identical across cohorts for columns chr/pos/ID.
4. **Strand orientation** (`03b_fix_strand.sh`, **disabled by design** — `FIX_STRAND="no"` in `config.sh`): an allele-based strand-flip step is included in the pipeline for completeness but turned off. Every SNP in every cohort classifies against each weight file's effect allele as ~86 % matched and ~14 % palindrome (A/T or C/G), with **no non-palindromic complement mismatches** — there is nothing for an allele-based strand fix to do on this dataset, and the only SNPs the script would touch are precisely the palindromes whose strand cannot be resolved from alleles alone (flipping them would invert their dosage encoding without justification). The step is kept in the repo, off, so it can be enabled if a future cohort or weight file shows real non-palindromic complement mismatches.
5. **Cross-cohort PCs** (`04b_compute_cross_cohort_pcs.sh`): 99,292 independent SNPs after excluding A/T and C/G palindromes + LD pruning (r² < 0.2 / 200 kb) on the pooled 18,898-sample merge.
6. **Cross-cohort z + raw canonical** (`07_assemble_output.R`): pooled (cross-cohort) mean/SD per trait → `PGI_trait_global_z`. The canonical bare column `PGI_trait` is set equal to `PGI_trait_raw` (raw PLINK2 `BETA_SUM`). With a single LOO weight file per trait across all cohorts, raw PGIs are already on a common scale, so no within-cohort standardization is applied — that would erase real cohort-level differences in genetic propensity (e.g. BASE-II sits ~0.2 SD higher on EA4 raw). `PGI_trait_global_z` is provided alongside for downstream code that wants z-units.

European-ancestry restriction is not part of this list: it happens at the merge
step, on the keep-lists step 04c produces (see the step table and
[`scripts/README_04c_ancestry.md`](scripts/README_04c_ancestry.md)).

---

## Leave-One-Out weight mapping (full matrix)

Single source of truth; any discrepancy between this table and `config.sh` is a bug.

| Cohort       | EA4                 | Cog   | NonCog | Height |
|--------------|---------------------|-------|--------|--------|
| BASE-II      | `EA4_excl_SHIP`     | `Cog` | `NCog` | `Height` |
| SHIP-0       | `EA4_excl_SHIP`     | `Cog` | `NCog` | `Height` |
| SHIP-Td      | `EA4_excl_SHIP`     | `Cog` | `NCog` | `Height` |
| SOEP-G       | `EA4_excl_SHIP`     | `Cog` | `NCog` | `Height` |
| TwinLife     | `EA4_excl_SHIP`     | `Cog` | `NCog` | `Height` |

Enforced by `get_ea4_trait_for_cohort()` and `get_noncog_trait_for_cohort()` in `config.sh`.

`NCog` is the published Malanchini NonCog. It is a different estimand from the
cognition-subtracted alternatives (Σβ +3.40) and shifts the cross-cohort raw
mean; that is a scaling property, not a validity problem, and it is neutralised
by within-cohort standardisation (`PGI_NonCog_within_z`, computed in
`03_Merge/merge.R`). Downstream, use the within-cohort z (or the global z) for
NonCog, never the raw cross-cohort mean.

---

## Verification: is this actually a harmonized PGI?

Run these on the cluster to confirm that all cohorts are scored on the identical SNP set.

```bash
cd "$TARDIS_PROJECT_ROOT"   # set by user_config.sh

# (A) Same number of SNPs across cohorts
for c in BASEII SHIP0 SHIPTD SOEP TWINLIFE; do
  wc -l output/geno/${c}_harmonized_shared.bim
done
# Expected: 6,878,689 five times

# (B) Identical chr/pos/ID columns
for c in BASEII SHIP0 SHIPTD SOEP TWINLIFE; do
  awk '{print $1"\t"$2"\t"$4}' output/geno/${c}_harmonized_shared.bim | md5sum
done
# Expected: one identical hash

# (C) Same NVAR scored per trait across cohorts
grep "variants processed" logs/score_*.log

# (D) Global harmonization statistics exist (step 07 ran)
cat output/diagnostics/global_pgi_statistics.tsv
```

(A)+(B) prove the SNP set is shared. (C) proves scoring used that shared set. (D) proves the cross-cohort standardization ran.

---

## Known data caveats

- **TwinLife's imputation panel (1000G Phase 3) differs** from the HRC v1.1 used by the other four cohorts. This is the reason TwinLife drives the shared-SNP intersection. It is managed — not removed — by the SNP-set restriction.
- **TwinLife EA4 contains one extreme value**, near −12.28 SD on `PGI_EA4`. It does not affect the reported bivariate correlations, but tail-sensitive analyses should winsorize or exclude it.
- **Within-family PGIs are not computed**: discovery summary statistics for a within-family EA GWAS are not available.

---

## Output

### Final per-cohort TSVs — `output/final/{COHORT}_PGI_PCs.tsv`

| Column(s)              | Description |
|------------------------|-------------|
| `IID`                  | Individual identifier (pseudonymous) |
| `PGI_{trait}_raw`      | Raw PLINK2 score (`BETA_SUM`) |
| `PGI_{trait}_global_z` | Cross-cohort z = `(raw − pooled_mean) / pooled_sd` |
| `PGI_{trait}`          | **Canonical PGI** — equal to `PGI_{trait}_raw`. See Cross-cohort harmonization §6 for why the bare column is raw, not z. |
| `PC1`–`PC20`           | Within-cohort ancestry PCs |
| `crossPC1`–`crossPC20` | Cross-cohort PCs (pooled) |

Traits: EA4, Cog, NonCog, Height → 4 traits × 3 PGI columns = 12 PGI columns.

> **Note on the downstream merged RDS.** `03_Merge/merge.R::normalize_pgi_columns()` repackages these columns at the merged-RDS layer: at that layer the bare `PGI_<trait>` is the analysis-friendly z (cross-cohort z in Harmonized mode), with `PGI_<trait>_raw` and `PGI_<trait>_within_z` preserved alongside. Analysis code in `04_GxE_Analysis/` reads the merged RDS and therefore sees z-units in the bare column; validation/sanity code should explicitly use `PGI_<trait>_raw`.

### Summary files
- `output/final/PGI_summary_statistics.tsv` — per-cohort descriptives for all PGIs
- `output/diagnostics/global_pgi_statistics.tsv` — global means/SDs used in step 07

---

## Source genotype data — on the SMB volumes

Each cohort's raw imputed genotype bfiles are sourced from a different
project share. Macs auto-mount these under `/Volumes/<volume>/`;
Workbench / ARC mounts them under `${SHARE_MOUNT_ROOT}/<volume>/` after
`biomount`. SMB URLs are `smb://<institution-share>/<volume>`.

| Cohort   | Share variable     | Path within volume                                  | bfile prefix                          |
|----------|--------------------|-----------------------------------------------------|---------------------------------------|
| BASE-II  | `${BASE2_SHARE}`   | `private/data/003_BASEII_geno/`                     | `LIFEBRAIN_BASEII_2015_AFFY`          |
| SHIP-0   | `${SHIP_SHARE}`    | `private/data/008_SHIP_processed_geno/`             | `Final-SHIP-0_R4a.chrall.dose`        |
| SHIP-Td  | `${SHIP_SHARE}`    | `private/data/008_SHIP_processed_geno/`             | `SHIP-Td-SHIP-Td_B2_merged.chrall`    |
| SOEP-G   | `${SOEP_SHARE}`    | `private/data/001_SOEP_processed_data/008_Genotype/`| `SOEP-G.b37.chrall`                   |
| TwinLife | `${TWINLIFE_SHARE}`| `private/data/Gendata/001_processed/`               | `TwinLife.b37.chrall`                 |

Each share variable is read from the environment; the defaults in the transfer
scripts are placeholders (`/path/to/<cohort>_share`).

Staging to the cluster is a one-time operation — see `transfer_geno_to_tardis.sh`
(macOS) or `transfer_geno_to_tardis_workbench.sh` (Workbench/ARC). On the
cluster they land in `data/geno/` of this project.

## On-cluster filesystem layout

The pipeline expects a self-contained project directory at
`${TARDIS_PROJECT_ROOT}` (set in `user_config.sh`) plus shared binaries
and reference data at `${TARDIS_BIN_DIR}` / `${TARDIS_LD_REF_DIR}`.
All data paths inside `config.sh` are anchored at `${PROJ_ROOT}`
(auto-detected from the config file's own location) so the project is
fully relocatable, and tool / LD-ref paths come from the env vars set
by `user_config.sh`.

The cluster is one beegfs filesystem mounted at `/mnt/beegfs` (~765 T);
the project root and the bin dir are different paths into the same pool,
so there is no reason to split data across them.

```
${TARDIS_PROJECT_ROOT}/                       # e.g. ${CLUSTER_PROJECT_ROOT}/80_years_.../genotype
├── config.sh                                 # paths, tools, trait/cohort mapping
├── user_config.sh                            # per-user overrides (gitignored)
├── _bootstrap_user_config.sh                 # interactive bootstrap helper
├── run_pipeline.sh                           # orchestrator (--from/--to/--step)
├── submit_post_sbayesr.slurm
├── README.md
├── scripts/                                  # 00–11 + setup_*.sh + transfer helpers
├── data/                                     # inputs — staged once; gitignored
│   ├── geno/                                 #   per-cohort imputed bfiles (5 prefixes)
│   ├── formatted_sumstats/                   #   per-trait, per-chr SBayesR-format
│   └── info_files/                           #   per-cohort imputation R² info files
├── output/                                   # run artifacts; gitignored
│   ├── sbayesr/                              #   per-chr SBayesR .snpRes / .parRes
│   ├── weights/                              #   scoring-format weights (chr:pos / A1 / BETA)
│   ├── geno/                                 #   {COHORT}_harmonized_shared.{bed,bim,fam}
│   ├── pcs/                                  #   within- and cross-cohort PC eigenvec/val
│   ├── ancestry/                             #   04c projected PCs + EUR keep-lists
│   ├── scores/                               #   PLINK2 .sscore files
│   ├── strand_fixes/                         #   03b flip lists (step disabled)
│   ├── diagnostics/                          #   global stats, allele freqs
│   └── final/                                #   per-cohort {COHORT}_PGI_PCs.tsv (deliverable)
└── logs/                                     # pipeline + SLURM logs (per run)

${TARDIS_BIN_DIR}/                            # e.g. ${CLUSTER_BIN_DIR}/bin
├── gctb                                      # SBayesR
├── plink                                     # 1.9
└── plink2                                    # 2.0 (harmonization, scoring, kinship)

${TARDIS_LD_REF_DIR}/                         # e.g. ${CLUSTER_BIN_DIR}/reference/ukb_50k_bigset_2.8M
└── ukb50k_shrunk_chr{1..22}_mafpt01.ldm.sparse[.bin]
```

`data/`, `output/`, and `logs/` are gitignored on the Mac side — they
live only on the machine that runs the pipeline.

## Mac ↔ cluster sync

| Direction | Tool | Notes |
|---|---|---|
| Mac → cluster (code only) | `bash transfer_pipeline_to_tardis.sh` | Pushes `01_Genotype/Harmonized/` from the monorepo to `${TARDIS_PROJECT_ROOT}`. Excludes `output/`, `logs/`, `data/`. `user_config.sh` IS synced (that's how config.sh on the cluster picks up your bin / LD-ref paths). |
| Mac/Workbench → cluster (genotypes + TwinLife R², one-time) | `bash transfer_geno_to_tardis.sh` *(macOS/SMB)* or `transfer_geno_to_tardis_workbench.sh` *(Workbench/ARC)* | Stages cohort bfiles into `data/geno/` and the TwinLife imputation-R² TSV into `data/info_files/TWINLIFE/`. Has retry logic and resumes via `transfer_geno.completed`. |
| Cluster → Mac/SMB (final outputs) | `rsync` direct (see [Phase 3](#phase-3--pull-results-back)) | Pulls `output/final/` and `output/diagnostics/` into the analysis project's `data/final/` mirror. The rest of `output/` is large and regenerable, so it is not mirrored. |

---

## Common pitfalls

- **SMB volumes unmount mid-transfer.** `transfer_geno_to_tardis.sh`
  detects this and tries to remount via the SMB URL list at the top of
  the script; if that fails it pauses and retries. If you see it loop
  on the same file, kill it (`Ctrl-C`), remount the share by hand
  (Finder → *Go* → *Connect to Server*), and re-run — progress resumes
  from `transfer_geno.completed`.
- **TwinLife R² source file** is a ~277 MiB pre-extracted aggregate TSV
  living on the TwinLife share at
  `private/data/Gendata/rawgendata/TwinLife_imputation_r2.tsv.gz`. The
  geno-transfer scripts pick it up from there automatically; it does not
  need re-extracting.
- **SBayesR can hit the 6 h SLURM wall.** The array submission
  (`scripts/01_run_sbayesr.slurm`) requests `--time=6:00:00`. EA4 LOO
  arrays usually finish in 2–4 h, but a slow node can push past 6 h;
  if any `sbr_*` job exits with `TIMEOUT`, resubmit just the missing
  arrays (the orchestrator's SKIP-if-complete check in
  `01_submit_all_sbayesr.sh` handles the rest).
- **SOEP-G GPG passphrase.** `setup_soep_geno.sh` requires a file with
  600 permissions (`chmod 600 ~/soep_passphrase.txt`). Delete the file
  after the one-time decryption + merge.
- **`config.sh` is location-anchored.** All paths inside it derive from
  `${PROJ_ROOT}` = the directory `config.sh` itself lives in. Don't
  symlink `config.sh` from another directory; the auto-detection will
  end up in the symlink source.
- **Passwordless SSH to the cluster** is required for the transfer scripts
  and the rsync in Phase 3. After running the user_config.sh setup, test
  with `ssh "$TARDIS_USER@$TARDIS_HOST" echo ok` — it must NOT prompt.
  If it does, run `ssh-copy-id "$TARDIS_USER@$TARDIS_HOST"` first.
- **`user_config.sh` is per-machine.** A copy lives on your Mac (loaded
  by every `transfer_*.sh` you run locally) AND on the cluster (rsync'd
  there by `transfer_pipeline_to_tardis.sh`; loaded by `config.sh`).
  Keep them aligned by always editing the Mac copy and re-running the
  pipeline transfer. Re-bootstrap any time by deleting the file and
  running any `transfer_*.sh` script — you'll be re-prompted.
- **Cohort-specific quirks** (SHIP-Td merge, SOEP encryption, TwinLife
  panel) are listed in [Cohort notes](#cohort-notes); a first-time
  runner should skim that section before starting Phase 1.
