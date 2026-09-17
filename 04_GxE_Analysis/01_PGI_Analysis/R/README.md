# R/

Shared helper functions for the 01_PGI_Analysis workflow, split by
responsibility. `GxE_Germany_analysis_helpers.R` is the single entry point:
each section script sources it, and it sources the `helpers_*.R` files below
in dependency order.

| File | Role |
|---|---|
| `helpers_formatting.R` | `z_scale`, `safe_numeric`, `z_within_group`, `cohort_bin` |
| `helpers_descriptives.R` | `kernel_z_score`, `kernel_moments`, `compute_attrition_flow` |
| `helpers_modeling.R` | `tidy_coefs`, `tidy_smooths`, `model_fit_row`, `compare_nested_models`, `lrt_pair`, `nested_lrt_table`, `assert_pgi_in_baseline`, `extract_*_pgi_effects`, `extract_pgi_slope`, `extract_lm_A1_curves`, `extract_a2_pgi_effects`, `diff_smooth_simul`, `diff_smooth_gender_simul`, `diff_smooth_pooled_gender_simul`, `compute_binned_slopes`, `within_cohort_slope_curves`, `focal_dispersion_check`, `cluster_vcov`, `tidy_coefs_cluster`, `focal_clust_vs_naive` |
| `helpers_plotting.R` | `plot_combined_east_west[_gender]`, `plot_gender_slopes_one_region`, `plot_diff_smooth`, `plot_cohort_coverage`, `plot_attrition_flow`, `plot_within_cohort_slopes` |
| `helpers_loo.R` | `cohort_groups()` (SHIP folds collapsed), `leave_one_cohort_out()` — the leave-one-study-out engine |
| `helpers_reporting.R` | `build_run_metadata`, `diff_smooth_significant_intervals`, `save_xlsx`, `multiblock_sheet`, `build_sheet_index`, `build_headline_table`, `build_gam_overview_table` (uses the `nonlinearity_conclusion` schema), `build_loo_table`, `write_section_xlsx`, `generate_section_overview_md`, `humanize_term/df/sheets` |
| `helpers_gender_loco_altpgi.R` | `rhs_without_re`, `loco_re_term`, `loco_pooled_coef`, `loco_pergender_gam`, `loco_pergender_curves`, `loco_b5_diffsmooth`, `altpgi_common_sample`, `altpgi_re_family`, `altpgi_slopes_family`, `altpgi_pergender_gam`, `re_pgi_slopes` |
| `helpers_descriptive_tables.R` | Descriptive-workbook building blocks: `fmt_msd`, `fmt_p`, `hedges_g`, `welch_p`, `group_cells`, `group_test`, `by_group3`, `desc_measure_row`, `build_desc_table[_by_group]`, `build_coverage_table`, `write_desc_table1/2`, `save_desc_workbook` |

`helpers_descriptive_tables.R` is not loaded by the wrapper; the
`00_Descriptives/` scripts source it directly.
