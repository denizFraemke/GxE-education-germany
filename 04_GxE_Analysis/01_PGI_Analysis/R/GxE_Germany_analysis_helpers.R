# R/GxE_Germany_analysis_helpers.R
#
# Thin wrapper: the helper function bodies live in sibling files in this
# same directory, ./helpers_*.R. This file is the single source-point each
# section's analysis/report script uses:
#   source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"),
#          local = TRUE)
#
# Scoping. The wrapper sources each helpers_*.R with local = TRUE,
# which evaluates the helper file in the env from which the helpers
# source() was called — i.e. the wrapper's own body frame. When a
# section script does
#   source(file.path(SCRIPT_DIR, "R", "GxE_Germany_analysis_helpers.R"),
#          local = TRUE)
# the wrapper's body frame is the main script's local frame, so the
# helper functions land there exactly as before. The dotted-prefixed
# temp vars below also live in that frame transiently and are rm()'d
# before the wrapper returns.

# ---- Resolve this wrapper's own location ----------------------------------
.wrapper_ofile <- NULL
for (.i in seq_len(sys.nframe())) {
  .frame <- sys.frame(.i)
  if (!is.null(.frame$ofile)) .wrapper_ofile <- .frame$ofile
}
if (is.null(.wrapper_ofile)) {
  stop("GxE_Germany_analysis_helpers.R wrapper: cannot resolve own ",
       "source-file location via sys.frames(). Source this file with ",
       "source(), not parse/eval.")
}
.wrapper_path <- normalizePath(.wrapper_ofile, mustWork = TRUE)
.helpers_root <- dirname(.wrapper_path)

# Sourced in dependency order: later files call functions defined earlier.
.helper_files <- c(
  "helpers_formatting.R",
  "helpers_descriptives.R",
  "helpers_modeling.R",
  "helpers_plotting.R",
  "helpers_loo.R",
  "helpers_reporting.R",
  "helpers_gender_loco_altpgi.R"
)
for (.f in .helper_files) {
  source(normalizePath(file.path(.helpers_root, .f), mustWork = TRUE),
         local = TRUE)
}

# ---- Clean up wrapper-internal variables ----------------------------------
rm(.wrapper_ofile, .wrapper_path, .helpers_root, .helper_files, .f, .i, .frame)

# ---- Sanity check ---------------------------------------------------------
# exists() walks the env chain, so this works regardless of whether the
# wrapper was sourced with local = TRUE (helpers land in caller's local
# frame) or local = FALSE (helpers land in global env).
stopifnot(
  exists("z_scale",                 mode = "function"),
  exists("safe_numeric",            mode = "function"),
  exists("kernel_z_score",          mode = "function"),
  exists("kernel_moments",          mode = "function"),
  exists("tidy_coefs",              mode = "function"),
  exists("tidy_smooths",            mode = "function"),
  exists("nested_lrt_table",        mode = "function"),
  exists("extract_pgi_slope",       mode = "function"),
  exists("diff_smooth_simul",       mode = "function"),
  exists("plot_combined_east_west", mode = "function"),
  exists("plot_diff_smooth",        mode = "function"),
  exists("plot_attrition_flow",     mode = "function"),
  exists("compute_attrition_flow",  mode = "function")
)
