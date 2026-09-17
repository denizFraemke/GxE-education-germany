# 00_setup/

Shared run-time setup, sourced by every section's analysis and report
script.

| File | Role |
|---|---|
| `constants.R` | `KERNEL_BW`, GAM `k` / basis, PGI outlier bound, descriptive cell floor, palettes |
| `load_data.R` | Locate and load the latest `Combined_Harmonized_*.rds` |
| `prepare_samples.R` | Build `df`, `dat_edu`, `dat_mob`, `dat_cluster`, `dat_cluster_mob`; define `by_mean` / `EDU_SD`; emit attrition |
| `run_context.R` | `get_run_ts()`, `data_root()`, `gxe_dir()`, `output_dir(section_id)` |
