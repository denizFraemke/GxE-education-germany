# paths.R — single source of truth for DATA_ROOT in the 03_Merge stage.
#
# Every script in 03_Merge resolves the analysis-share root through
# data_root(): set the DATA_ROOT environment variable to override it
# (e.g. point at a replication folder); otherwise the default below is
# used.
#
# Change the default in THIS one place only — never inline it in
# individual scripts (that is how a path edit gets missed somewhere and
# a run silently reads the wrong data).
data_root <- function() {
  Sys.getenv(
    "DATA_ROOT",
    unset = "/path/to/data_root"
  )
}
