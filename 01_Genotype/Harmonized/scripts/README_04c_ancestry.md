# 04c — European-ancestry keep-lists (1000G projection)

Reference-based ancestry exclusion applied at the merge step. The merge
(`03_Merge`) consumes the per-cohort keep-lists via `filter_to_eur()`; this step
produces them.

## Method

1. **Reference** — 1000 Genomes phase 3, GRCh37/b37 (autosomes), the build that
   matches our harmonized cohort genotypes. Public download (see below).
2. **`04c_ancestry_1kg.sh`** — builds the reference PLINK set (biallelic ACGT
   SNPs, `chr:pos` IDs, MAF filter, long-range-LD excluded, strand-ambiguous
   SNPs dropped, LD-pruned), computes reference PCA with allele weights, then
   projects the reference samples and each cohort onto those PCs
   (`--score … variance-standardize`, identical projection scale for all).
3. **`04c_classify_eur.R`** — in the first *K* projected PCs, computes each
   individual's Mahalanobis distance to every 1000G super-population centroid
   (each with that population's own covariance) and **keeps individuals whose
   nearest centroid is EUR and whose distance to the EUR centroid is within the
   chi-square cutoff (df = K)**. Defaults: `K = 6`, cutoff = χ²₀.₉₉₉. Both are
   parameters and are echoed to the log.

Output: `EUR_keep_IDs_<KEY>.tsv` (single `IID` column) per cohort, to be placed
in `data/1Kg_exclusion_IDs/`. Keys: `BASEII`, `SHIP0`, `SHIPTd`, `TwinLife`.

## Per-cohort provenance

- **BASE-II, SHIP (SHIP-0 + SHIP-Td), TwinLife** — classified in-house by the
  1000G projection above.
- **SOEP** — the SOEP genotypes are not available to project, so
  `EUR_keep_IDs_SOEP.tsv` comes from a European-ancestry classification carried
  out on the strict-QC (chrall) SOEP release. The scored sample here is the
  mild-QC (ext_chrall) release, so the merge keeps the on-list Europeans **plus**
  the ext-only individuals absent from chrall (unassessable, ~98% European per
  the cohort profile) and drops only the chrall individuals classified
  non-European (`geno_qc` tags each row `strict` / `mild_only`).

## Download the reference (public, run on the cluster)

```bash
REF="$HOME/reference/1000G_phase3_b37"      # set your path
mkdir -p "$REF" && cd "$REF"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
wget -c "${BASE}/integrated_call_samples_v3.20130502.ALL.panel"
for chr in $(seq 1 22); do
  wget -c "${BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  wget -c "${BASE}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
done
```

## Run

```bash
# on the cluster, after the reference download + step 03 (harmonized bfiles):
REF_1KG_VCF_DIR="$HOME/reference/1000G_phase3_b37" bash scripts/04c_ancestry_1kg.sh
EXISTING_KEEP_DIR=<path to current keep-lists> \
  Rscript scripts/04c_classify_eur.R          # EXISTING_KEEP_DIR is optional (concordance check)
```

## Validation

The R script prints, per cohort, the EUR-kept / dropped counts and — if
`EXISTING_KEEP_DIR` is given — the agreement against the keep-lists in that
directory. Confirm high concordance (and inspect any large new-only /
dropped-vs-old sets) before regenerating the merge. Tune `N_PC` / `EUR_Q` only
with a stated rationale.
